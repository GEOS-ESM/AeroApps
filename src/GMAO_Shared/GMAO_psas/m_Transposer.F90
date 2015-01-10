!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_Transposer - Transpose a distributed vector to another
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_Transposer
      implicit none
      private	! except

      public :: Transposer		! The class data structure

      public :: init ,Transposer_init	! intialize a Transposer
      public :: clean,Transposer_clean	! clean a Transposer

      public :: new	! create an object as a target
      public :: delete	! delete an object as a pointer

      public :: sendSize		! "send" size
      public :: recvSize		! "recv" size
      public :: getRange_send		! get-range at the "send" side
      public :: getRange_recv		! get-range at the "recv" side

      public :: ptr_sendCounts	! reference of %scounts
      public :: ptr_sendDispls	! reference of %sdispls
      public :: ptr_recvCounts	! reference of %rcounts
      public :: ptr_recvDispls	! reference of %rdispls

    type Transposer
      private
      integer :: myPE
      integer,pointer,dimension(:) :: scounts
      integer,pointer,dimension(:) :: sdispls
      integer,pointer,dimension(:) :: rcounts
      integer,pointer,dimension(:) :: rdispls
    end type Transposer

    interface Transposer_init ; module procedure init_ ; end interface
    interface init ; module procedure init_ ; end interface
    interface Transposer_clean; module procedure clean_; end interface
    interface clean; module procedure clean_; end interface

    interface new   ; module procedure new_   ; end interface
    interface delete; module procedure delete_; end interface

    interface sendSize; module procedure ssize_; end interface
    interface recvSize; module procedure rsize_; end interface

    interface getRange_send; module procedure srange_; end interface
    interface getRange_recv; module procedure rrange_; end interface

    interface ptr_sendCounts; module procedure	&
	ptr_scounts_; end interface
    interface ptr_sendDispls; module procedure	&
	ptr_sdispls_; end interface
    interface ptr_recvCounts; module procedure	&
	ptr_rcounts_; end interface
    interface ptr_recvDispls; module procedure	&
	ptr_rdispls_; end interface

! !REVISION HISTORY:
! 	24Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_Transposer'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize a transposer
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(trans,scounts,comm,stat)
      use m_mpif90,only : MP_comm_size,MP_comm_rank
      use m_mpif90,only : MP_perr
      use m_mpif90,only : MP_type
      use m_die,only : perr,die
      use m_mall,only : mall_ison,mall_mci
      implicit none
      type(Transposer),intent(out) :: trans
      integer,dimension(0:),intent(in) :: scounts
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	24Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: nPE,myPE
  integer :: mesgType
  integer :: i,ier

  call MP_comm_size(comm,nPE,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MP_comm_size()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  if(nPE/=size(scounts)) then
    call perr(myname_,'nPE',nPE,'size(scounts)',size(scounts))
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

  allocate( trans%scounts(0:nPE-1),trans%sdispls(0:nPE-1),	&
	    trans%rcounts(0:nPE-1),trans%rdispls(0:nPE-1),	&
	    stat=ier	)
	if(ier/=0) then
	  call MP_perr(myname_,'allocate()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

	if(mall_ison()) then
	  call mall_mci(trans%scounts,myname)
	  call mall_mci(trans%sdispls,myname)
	  call mall_mci(trans%rcounts,myname)
	  call mall_mci(trans%rdispls,myname)
	endif

  call MP_comm_rank(comm,myPE,ier)
	if(ier/=0) then
	  call MP_perr(myname_,'MP_comm_rank()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
	
  trans%myPE   =myPE
  trans%scounts=scounts

  mesgType=MP_type(trans%scounts(0))
  call MPI_alltoall(trans%scounts,1,mesgType,   &
                    trans%rcounts,1,mesgType,comm,ier)
        if(ier/=0) then
          call MP_perr(myname_,'MPI_alltoall()',ier)
          if(.not.present(stat)) call die(myname_)
          stat=ier
          return
        endif

  trans%sdispls(0)=0
  do i=1,nPE-1
    trans%sdispls(i)=trans%sdispls(i-1)+trans%scounts(i-1)
  end do

  trans%rdispls(0)=0
  do i=1,nPE-1
    trans%rdispls(i)=trans%rdispls(i-1)+trans%rcounts(i-1)
  end do
end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean a Transposer
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(trans,stat)
      use m_die,only : die
      use m_mall,only : mall_ison,mall_mco
      implicit none
      type(Transposer),intent(inout) :: trans
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	24Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

  if(present(stat)) stat=0

	if(mall_ison()) then
	  call mall_mco(trans%scounts,myname)
	  call mall_mco(trans%sdispls,myname)
	  call mall_mco(trans%rcounts,myname)
	  call mall_mco(trans%rdispls,myname)
	endif

  deallocate( trans%scounts,trans%sdispls,	&
	      trans%rcounts,trans%rdispls, stat=ier)
	if(ier/=0) then
	  if(.not.present(stat)) call die(myname_,'deallocate()',ier)
	  stat=ier
	  return
	endif

end subroutine clean_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_scounts_ - get the pointer to %scounts
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_scounts_(trans)
      implicit none
      type(Transposer),intent(in) :: trans
      integer,pointer,dimension(:) :: ptr_scounts_

! !REVISION HISTORY:
! 	24Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_scounts_'

  ptr_scounts_ => trans%scounts

end function ptr_scounts_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_sdispls_ - get the pointer to %sdispls
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_sdispls_(trans)
      implicit none
      type(Transposer),intent(in) :: trans
      integer,pointer,dimension(:) :: ptr_sdispls_

! !REVISION HISTORY:
! 	24Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_sdispls_'

  ptr_sdispls_ => trans%sdispls

end function ptr_sdispls_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_rcounts_ - get the pointer to %rcounts
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_rcounts_(trans)
      implicit none
      type(Transposer),intent(in) :: trans
      integer,pointer,dimension(:) :: ptr_rcounts_

! !REVISION HISTORY:
! 	24Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_rcounts_'

  ptr_rcounts_ => trans%rcounts

end function ptr_rcounts_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_rdispls_ - get the pointer to %rdispls
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_rdispls_(trans)
      implicit none
      type(Transposer),intent(in) :: trans
      integer,pointer,dimension(:) :: ptr_rdispls_

! !REVISION HISTORY:
! 	24Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_rdispls_'

  ptr_rdispls_ => trans%rdispls

end function ptr_rdispls_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ssize_ - the size at the send side
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ssize_(trans)
      implicit none
      type(Transposer),intent(in) :: trans
      integer :: ssize_

! !REVISION HISTORY:
! 	31Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ssize_'
  integer :: mPE

  mPE=size(trans%scounts)-1
  ssize_=trans%sdispls(mPE)+trans%scounts(mPE)

end function ssize_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rsize_ - the size at the recv side
!
! !DESCRIPTION:
!
! !INTERFACE:

    function rsize_(trans)
      implicit none
      type(Transposer),intent(in) :: trans
      integer :: rsize_

! !REVISION HISTORY:
! 	31Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rsize_'
  integer :: mPE

  mPE=size(trans%rcounts)-1
  rsize_=trans%rdispls(mPE)+trans%rcounts(mPE)

end function rsize_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: srange_ - get the range at the "send" side
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine srange_(trans,lbound,ubound,iPE)
      implicit none
      type(Transposer),intent(in) :: trans
      integer,intent(out) :: lbound,ubound
      integer,optional,intent(in) :: iPE

! !REVISION HISTORY:
! 	27Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::srange_'
  integer :: iPE_

  iPE_=trans%myPE
  if(present(iPE)) iPE_=iPE

  lbound=trans%sdispls(iPE_)+1
  ubound=trans%sdispls(iPE_)+trans%scounts(iPE_)

end subroutine srange_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rrange_ - get the range at the "recv" side
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine rrange_(trans,lbound,ubound,iPE)
      implicit none
      type(Transposer),intent(in) :: trans
      integer,intent(out) :: lbound,ubound
      integer,optional,intent(in) :: iPE

! !REVISION HISTORY:
! 	27Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rrange_'
  integer :: iPE_

  iPE_=trans%myPE
  if(present(iPE)) iPE_=iPE

  lbound=trans%rdispls(iPE_)+1
  ubound=trans%rdispls(iPE_)+trans%rcounts(iPE_)

end subroutine rrange_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: new_ - create an object as a target
!
! !DESCRIPTION:
!
! !INTERFACE:

    function new_(mold,stat)
      use m_die, only : die
      use m_mall,only : mall_ison,mall_ci
      implicit none
      type(Transposer),pointer :: mold
      integer,optional,intent(out) :: stat
      type(Transposer),pointer :: new_

! !REVISION HISTORY:
! 	28Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::new_'
  integer :: ier
  type(Transposer),pointer :: trans

  if(present(stat)) stat=0
  allocate(trans,stat=ier)
	if(ier/=0) then
	  if(.not.present(stat)) call die(myname_,'allocate()',ier)
	  stat=ier
	  return
	endif

	if(mall_ison()) call mall_ci(1,myname)

  new_ => trans
  nullify(trans)	! to prevent the compiler touching the memory.
end function new_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: delete_ - delete an object as a pointer
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine delete_(tran,stat)
      use m_die, only : die
      use m_mall,only : mall_ison,mall_co
      implicit none
      type(Transposer),pointer :: tran
      integer,optional,intent(out) :: stat


! !REVISION HISTORY:
! 	28Mar00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::delete_'
  integer :: ier

  if(present(stat)) stat=0

	if(mall_ison()) call mall_co(1,myname)

  deallocate(tran,stat=ier)
	if(ier/=0) then
	  if(.not.present(stat)) call die(myname_,'deallocate()',ier)
	  stat=ier
	  return
	endif

end subroutine delete_

end module m_Transposer
