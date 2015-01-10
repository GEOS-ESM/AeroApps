!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_Partition - Partition-triple (indx,tran,coll)
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_Partition
      use m_Collector ,only : Collector
      use m_Transposer,only : Transposer
      implicit none
      private	! except

      public :: Partition		! Class data structure
      public :: Partition_init ,init	! initialize
      public :: Partition_clean,clean	! clean
      public :: localSize		! local size
      public :: globalSize		! global size
      public :: distrSize		! distributed size
      public :: ptr_indx		! refering to %indx
      public :: ptr_coll		! refering to %coll
      public :: ptr_tran		! refering to %tran
      
      type Partition
	private
        integer,pointer,dimension(:) :: indx ! pre-distribution sorting
        type(Collector) ,pointer :: coll
        type(Transposer),pointer :: tran
      end type Partition

      interface Partition_init; module procedure	&
	init_; end interface
      interface init; module procedure init_; end interface
        
      interface Partition_clean; module procedure	&
	clean_; end interface
      interface clean; module procedure clean_; end interface

      interface  localSize; module procedure lsize_; end interface
      interface globalSize; module procedure gsize_; end interface
      interface  distrSize; module procedure dsize_; end interface

      interface ptr_indx; module procedure ptr_indx_; end interface
      interface ptr_coll; module procedure ptr_coll_; end interface
      interface ptr_tran; module procedure ptr_tran_; end interface
        
! !REVISION HISTORY:
! 	09Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_Partition'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - Create a Partition-triple for a given vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(part,attr,ptnr,comm)
      use m_die,only : die
      use m_mall,only : mall_ison,mall_ci,mall_mci
      use m_Attributes,only : Attributes
      use m_Attributes,only : lsize

      use m_Partitioner,only : Partitioner
      use m_Partitioner,only : Partitioner_partition

      implicit none
      type(Partition)   ,intent(out) :: part
      type(Attributes)  ,intent(in)  :: attr  ! to be partitioned
      type(Partitioner) ,intent(in)  :: ptnr
      integer           ,intent(in)  :: comm  ! Communicator

! !REVISION HISTORY:
! 	09Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: lsize_
  integer :: ier

  lsize_=lsize(attr)
  allocate(part%indx(lsize_),part%tran,part%coll,stat=ier)
        if(ier/=0) call die(myname_,'allocate()',ier)
        if(mall_ison()) then
          call mall_mci(part%indx,myname)
          call mall_ci(2,myname)
        endif

  call Partitioner_partition(attr,      &
        part%indx,part%tran,part%coll, ptnr,comm)

end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - Clean a Partition-triple
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(part)
      use m_die,only : die
      use m_mall,only : mall_ison,mall_co,mall_mco
      use m_Collector ,only : Collector_clean
      use m_Transposer,only : Transposer_clean

      implicit none
      type(Partition)   ,intent(inout) :: part

! !REVISION HISTORY:
! 	09Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier
  
  call Collector_clean(part%coll)
  call Transposer_clean(part%tran)

        if(mall_ison()) then
          call mall_co(2,myname)
          call mall_mco(part%indx,myname)
        endif
        
  deallocate(part%coll,part%tran,part%indx,stat=ier)
        if(ier/=0) call die(myname_,'deallocate()',ier)
end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: lsize_ - local size of a partition
!
! !DESCRIPTION:
!
! !INTERFACE:

    function lsize_(part)
      use m_Collector,only : localSize
      implicit none
      type(Partition),intent(in) :: part
      integer :: lsize_

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::lsize_'

  lsize_ = localSize(part%coll)

end function lsize_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gsize_ - global size of a partition
!
! !DESCRIPTION:
!
! !INTERFACE:

    function gsize_(part)
      use m_Collector,only : globalSize
      implicit none
      type(Partition),intent(in) :: part
      integer :: gsize_

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gsize_'

  gsize_=globalSize(part%coll)

end function gsize_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: dsize_ - distributed size of a partition
!
! !DESCRIPTION:
!
! !INTERFACE:

    function dsize_(part)
      use m_Transposer,only : recvSize
      implicit none
      type(Partition),intent(in) :: part
      integer :: dsize_

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::dsize_'

  dsize_ = recvSize(part%tran)

end function dsize_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_indx_ - refering to %indx
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_indx_(part)
      implicit none
      type(Partition),intent(in) :: part
      integer,pointer,dimension(:) :: ptr_indx_

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_indx_'

  ptr_indx_ => part%indx

end function ptr_indx_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_tran_ - refering to %tran
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_tran_(part)
      use m_Transposer,only : Transposer
      implicit none
      type(Partition),intent(in) :: part
      type(Transposer),pointer :: ptr_tran_

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_tran_'

  ptr_tran_ => part%tran

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

    function ptr_coll_(part)
      use m_Collector,only : Collector
      implicit none
      type(Partition),intent(in) :: part
      type(Collector),pointer :: ptr_coll_

! !REVISION HISTORY:
! 	11Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_coll_'

  ptr_coll_ => part%coll

end function ptr_coll_

end module m_Partition
