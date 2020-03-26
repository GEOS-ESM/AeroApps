!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_Navigator - Array of pointers
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_Navigator
      implicit none
      private	! except

      public :: Navigator		! The class data structure
      public :: Navigator_init!,init	! initialize an object
      public :: clean			! clean an object
      public :: new			! create an object
      public :: delete			! delete an object
      public :: lsize			! the true size
      public :: msize			! the maximum size
      public :: resize			! adjust the true size
      public :: get			! get an entry
      public :: inquire			! inquire memory information

      public :: define			! define a navigator by counts
      public :: weighted_define		! define a navigator by weights

      public :: allgather		! allgather navigators

      public :: ptr_displs		! referencing %displs(:)
      public :: ptr_counts		! referencing %counts(:)

    type Navigator
      private
      integer :: lsize	! true size, if over-dimensioned
      integer,pointer,dimension(:) :: displs
      integer,pointer,dimension(:) :: counts
    end type Navigator

    interface Navigator_init; module procedure	&
	init_,		&
	initbykeys_,	&
	initby2keys_,	&
	initby3keys_,	&
	initby4keys_,	&
	initby5keys_; end interface
    interface init  ; module procedure		&
	init_,		&
	initbykeys_,	&
	initby2keys_,	&
	initby3keys_,	&
	initby4keys_,	&
	initby5keys_; end interface

    interface Navigator_define  ; module procedure	&
	define_  ; end interface
    interface define  ; module procedure define_  ; end interface
    interface clean ; module procedure clean_ ; end interface

    interface lsize ; module procedure lsize_ ; end interface
    interface msize ; module procedure msize_ ; end interface
    interface resize; module procedure resize_; end interface
    interface get   ; module procedure get_   ; end interface

    interface allgather; module procedure allgather_; end interface

    interface ptr_displs; module procedure	&
	ptr_displs_; end interface
    interface ptr_counts; module procedure	&
	ptr_counts_; end interface

    interface    new; module procedure    new_; end interface
    interface delete; module procedure delete_; end interface

    interface weighted_define; module procedure	&
	weighted_define_; end interface

    interface inquire; module procedure inquire_; end interface

! !REVISION HISTORY:
! 	22May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_Navigator'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize a navigator
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(nav,lsize,stat)
      use m_mall,only : mall_ison,mall_mci
      use m_die ,only : die,perr
      implicit none
      type(Navigator),intent(out) :: nav	! the object
      integer,intent(in) :: lsize		! nominal size
      integer,optional,intent(out) :: stat	! status

! !REVISION HISTORY:
! 	22May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: ier

  if(present(stat)) stat=0

  allocate(nav%displs(lsize),nav%counts(lsize),stat=ier)
	if(ier/=0) then
	  call perr(myname_,'allocate()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
	if(mall_ison()) then
	  call mall_mci(nav%displs,myname)
	  call mall_mci(nav%counts,myname)
	endif

  nav%lsize=lsize
end subroutine init_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initbykeys_ - define navigator with sorted keys
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine initbykeys_(nav,keys,stat)
      use m_die ,only : die,perr
      use m_mall,only : mall_ison,mall_mci,mall_mco
      implicit none
      type(Navigator),intent(out) :: nav
      integer,dimension(:),intent(in) :: keys
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	23Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initbykeys_'
  integer,allocatable,dimension(:) :: counts
  integer :: ier
  integer :: i,n
  integer :: nks
  integer :: key_i

  if(present(stat)) stat=0

  n=size(keys)

  if(n>0) then

	allocate(counts(n),stat=ier)
		if(ier/=0) then
		  call perr(myname_,'allocate()',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
		if(mall_ison()) call mall_mci(counts,myname)

    nks=1
    key_i=keys(1)	! save the key

    counts(nks)=1

    do i=2,n
      if(keys(i)/=key_i) then
	nks=nks+1	! next sounding
	counts(nks)=0	! initialized the new counter
	key_i=keys(i)	! save the key
      endif
      counts(nks)=counts(nks)+1
    end do

    call init_(nav,nks,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'init_()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

    call define_(nav,counts=counts(1:nks))

		if(mall_ison()) call mall_mco(counts,myname)
	deallocate(counts,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'deallocate()',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif

  else

    call init_(nav,0,stat=ier)

	if(ier/=0) then
	  call perr(myname_,'init_()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
  endif

end subroutine initbykeys_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initby2keys_ - define navigator with 2 sorted keys
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine initby2keys_(nav,key1,key2,stat)
      use m_die ,only : die,perr
      use m_mall,only : mall_ison,mall_mci,mall_mco
      implicit none
      type(Navigator),intent(out) :: nav
      integer,dimension(:),intent(in) :: key1
      integer,dimension(:),intent(in) :: key2
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
!       03Jul00 - lpc
!               - 'key1_i=key2(i)  ! save key2' is changed into
!                 'key2_i=key2(i)  ! save key2'
!               - after nks is determined, use counts(1:nks) 
!                 to define a new navigator 'nav'
! 	23Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initby2keys_'
  integer,allocatable,dimension(:) :: counts
  integer :: ier
  integer :: i,n
  integer :: nks
  integer :: key1_i,key2_i

  if(present(stat)) stat=0

  n=size(key1)
  if(n>size(key2)) call die(myname_,'invalid size(key2)',n)

  if(n>0) then

	allocate(counts(n),stat=ier)
		if(ier/=0) then
		  call perr(myname_,'allocate(counts)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
		if(mall_ison()) call mall_mci(counts,myname)

    nks=1
    key1_i=key1(1)	! save the key
    key2_i=key2(1)

    counts(nks)=1

    do i=2,n
      if(key1(i)/=key1_i.or.key2(i)/=key2_i) then
	nks=nks+1	! next sounding
	counts(nks)=0	! initialized the new counter
	key1_i=key1(i)	! save key1
	key2_i=key2(i)	! save key2; by lpc 0703,2000
      endif
      counts(nks)=counts(nks)+1
    end do

    call init_(nav,nks,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'init_()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif


    call define_(nav,counts=counts(1:nks))

		if(mall_ison()) call mall_mco(counts,myname)
	deallocate(counts,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'deallocate(counts)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif

  else

    call init_(nav,0,stat=ier)

	if(ier/=0) then
	  call perr(myname_,'init_()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
  endif

end subroutine initby2keys_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initby3keys_ - define navigator with 3 sorted keys
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine initby3keys_(nav,key1,key2,key3,stat)
      use m_die ,only : die,perr
      use m_mall,only : mall_ison,mall_mci,mall_mco
      implicit none
      type(Navigator),intent(out) :: nav
      integer,dimension(:),intent(in) :: key1
      integer,dimension(:),intent(in) :: key2
      integer,dimension(:),intent(in) :: key3
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
!	11Sep00	- Jing Guo
!		. Extended from 2 keys to 3 keys.
!       03Jul00 - lpc
!               - 'key1_i=key2(i)  ! save key2' is changed into
!                 'key2_i=key2(i)  ! save key2'
!               - after nks is determined, use counts(1:nks) 
!                 to define a new navigator 'nav'
! 	23Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initby3keys_'
  integer,allocatable,dimension(:) :: counts
  integer :: ier
  integer :: i,n
  integer :: nks
  integer :: key1_i,key2_i,key3_i

  if(present(stat)) stat=0

  n=size(key1)
  if(n>size(key2)) call die(myname_,'invalid size(key2)',n)
  if(n>size(key3)) call die(myname_,'invalid size(key3)',n)

  if(n>0) then

	allocate(counts(n),stat=ier)
		if(ier/=0) then
		  call perr(myname_,'allocate(counts)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
		if(mall_ison()) call mall_mci(counts,myname)

    nks=1
    key1_i=key1(1)	! save the key
    key2_i=key2(1)
    key3_i=key3(1)

    counts(nks)=1

    do i=2,n
      if( key1(i)/=key1_i .or.	&
	  key2(i)/=key2_i .or.	&
	  key3(i)/=key3_i ) then
	nks=nks+1	! next sounding
	counts(nks)=0	! initialized the new counter
	key1_i=key1(i)	! save key1
	key2_i=key2(i)	! save key2; by lpc 0703,2000
	key3_i=key3(i)	! save key3
      endif
      counts(nks)=counts(nks)+1
    end do

    call init_(nav,nks,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'init_()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

    call define_(nav,counts=counts(1:nks))

		if(mall_ison()) call mall_mco(counts,myname)
	deallocate(counts,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'deallocate(counts)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif

  else

    call init_(nav,0,stat=ier)

	if(ier/=0) then
	  call perr(myname_,'init_()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
  endif

end subroutine initby3keys_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initby4keys_ - define navigator with 4 sorted keys
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine initby4keys_(nav,key1,key2,key3,key4,stat)
      use m_die ,only : die,perr
      use m_mall,only : mall_ison,mall_mci,mall_mco
      implicit none
      type(Navigator),intent(out) :: nav
      integer,dimension(:),intent(in) :: key1
      integer,dimension(:),intent(in) :: key2
      integer,dimension(:),intent(in) :: key3
      integer,dimension(:),intent(in) :: key4
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
!	11Sep00	- Jing Guo
!		. Extended from 3 keys to 4 keys.
!       03Jul00 - lpc
!               - 'key1_i=key2(i)  ! save key2' is changed into
!                 'key2_i=key2(i)  ! save key2'
!               - after nks is determined, use counts(1:nks) 
!                 to define a new navigator 'nav'
! 	23Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initby4keys_'
  integer,allocatable,dimension(:) :: counts
  integer :: ier
  integer :: i,n
  integer :: nks
  integer :: key1_i,key2_i,key3_i,key4_i

  if(present(stat)) stat=0

  n=size(key1)
  if(n>size(key2)) call die(myname_,'invalid size(key2)',n)
  if(n>size(key3)) call die(myname_,'invalid size(key3)',n)
  if(n>size(key4)) call die(myname_,'invalid size(key4)',n)

  if(n>0) then

	allocate(counts(n),stat=ier)
		if(ier/=0) then
		  call perr(myname_,'allocate(counts)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
		if(mall_ison()) call mall_mci(counts,myname)

    nks=1
    key1_i=key1(1)	! save the key
    key2_i=key2(1)
    key3_i=key3(1)
    key4_i=key4(1)

    counts(nks)=1

    do i=2,n
      if( key1(i)/=key1_i .or.	&
	  key2(i)/=key2_i .or.	&
	  key3(i)/=key3_i .or.	&
	  key4(i)/=key4_i ) then
	nks=nks+1	! next sounding
	counts(nks)=0	! initialized the new counter
	key1_i=key1(i)	! save key1
	key2_i=key2(i)	! save key2; by lpc 0703,2000
	key3_i=key3(i)	! save key3
	key4_i=key4(i)	! save key4
      endif
      counts(nks)=counts(nks)+1
    end do

    call init_(nav,nks,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'init_()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

    call define_(nav,counts=counts(1:nks))

		if(mall_ison()) call mall_mco(counts,myname)
	deallocate(counts,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'deallocate(counts)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif

  else

    call init_(nav,0,stat=ier)

	if(ier/=0) then
	  call perr(myname_,'init_()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
  endif

end subroutine initby4keys_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initby5keys_ - define navigator with 4 sorted keys
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine initby5keys_(nav,key1,key2,key3,key4,key5,stat)
      use m_die ,only : die,perr
      implicit none
      type(Navigator),intent(out) :: nav
      integer,dimension(:),intent(in) :: key1
      integer,dimension(:),intent(in) :: key2
      integer,dimension(:),intent(in) :: key3
      integer,dimension(:),intent(in) :: key4
      integer,dimension(:),intent(in) :: key5
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
!	11Sep00	- Jing Guo
!		. Extended from 3 keys to 4 keys.
!       03Jul00 - lpc
!               - 'key1_i=key2(i)  ! save key2' is changed into
!                 'key2_i=key2(i)  ! save key2'
!               - after nks is determined, use counts(1:nks) 
!                 to define a new navigator 'nav'
! 	23Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initby5keys_'
  integer,allocatable,dimension(:) :: counts
  integer :: ier
  integer :: i,n
  integer :: nks
  integer :: key1_i,key2_i,key3_i,key4_i,key5_i

  if(present(stat)) stat=0

  n=size(key1)
  if(n>size(key2)) call die(myname_,'invalid size(key2)',n)
  if(n>size(key3)) call die(myname_,'invalid size(key3)',n)
  if(n>size(key4)) call die(myname_,'invalid size(key4)',n)
  if(n>size(key5)) call die(myname_,'invalid size(key5)',n)

  if(n>0) then

	allocate(counts(n),stat=ier)
		if(ier/=0) then
		  call perr(myname_,'allocate(counts)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif

    nks=1
    key1_i=key1(1)	! save the key
    key2_i=key2(1)
    key3_i=key3(1)
    key4_i=key4(1)
    key5_i=key5(1)

    counts(nks)=1

    do i=2,n
      if( key1(i)/=key1_i .or.	&
	  key2(i)/=key2_i .or.	&
	  key3(i)/=key3_i .or.	&
	  key4(i)/=key4_i .or.	&
	  key5(i)/=key5_i ) then
	nks=nks+1	! next sounding
	counts(nks)=0	! initialized the new counter
	key1_i=key1(i)	! save key1
	key2_i=key2(i)	! save key2; by lpc 0703,2000
	key3_i=key3(i)	! save key3
	key4_i=key4(i)	! save key4
	key5_i=key5(i)	! save key5
      endif
      counts(nks)=counts(nks)+1
    end do

    call init_(nav,nks,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'init_()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

    call define_(nav,counts=counts(1:nks))

	deallocate(counts,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'deallocate(counts)',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif

  else

    call init_(nav,0,stat=ier)

	if(ier/=0) then
	  call perr(myname_,'init_()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
  endif

end subroutine initby5keys_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: define_ - define navigator components
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine define_(nav,counts,displs)
      use m_die,only : die
      implicit none
      type(Navigator),intent(inout) :: nav	! the object
      integer,optional,dimension(:),intent(in) :: counts
      integer,optional,dimension(:),intent(in) :: displs

! !REVISION HISTORY:
! 	15Jun00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::define_'
  integer :: k,l,m,n

! Case 1:  See also, Case 2
!
!     call Navigator_init(nav,n)
!	counts=>ptr_counts(nav)
!	counts(:)=...; nullify(counts)
!     call Navigator_define(nav)
!
! Case 2:  Both Case 1 and Case 2 allow safe definition of %displs
!
!	lns(:)=...
!     call Navigator_init(nav,n)
!     call Navigator_define(nav,counts=lns)
!
! Case 3:  It allows non-sequential displs(:) definitions
!
!	lns(:)=...
!	lcs(:)=..
!     call Navigator_init(nav,n)
!     call Navigator_define(nav,counts=lns,displs=lcs(:)-1)

  m=msize_(nav)

  if(present(counts)) then
    l=size(counts)
    if(l>m) call die(myname_,'size(counts)',l)

    nav%counts(1:l)=counts
    call resize_(nav,l)
  endif

  if(present(displs)) then
    l=size(displs)
    n=lsize_(nav)
    if(l>n) call die(myname_,'size(displs)',l)

    nav%displs(1:l)=displs

  else
	! Define displs(:) with a default sequential definition.

    n=lsize_(nav)
    if(n<=0) return

    nav%displs(1)=0
    do k=2,n
      nav%displs(k)=nav%displs(k-1)+nav%counts(k-1)
    end do
  endif

end subroutine define_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean a navigator
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(nav,stat)
      use m_mall,only : mall_ison,mall_mco
      use m_die ,only : die,perr
      implicit none
      type(Navigator),intent(inout) :: nav	! the object
      integer,optional,intent(out) :: stat	! status

! !REVISION HISTORY:
! 	22May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

  if(present(stat)) stat=0

	if(mall_ison()) then
	  call mall_mco(nav%displs,myname)
	  call mall_mco(nav%counts,myname)
	endif
  deallocate(nav%displs,nav%counts,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'deallocate()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  nav%lsize=0
end subroutine clean_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: lsize_ - return the true size
!
! !DESCRIPTION:
!
! !INTERFACE:

    function lsize_(nav)
      implicit none
      type(Navigator),intent(in) :: nav
      integer :: lsize_

! !REVISION HISTORY:
! 	22May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::lsize_'

  lsize_=nav%lsize
end function lsize_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: msize_ - return the maximum size
!
! !DESCRIPTION:
!
! !INTERFACE:

    function msize_(nav)
      implicit none
      type(Navigator),intent(in) :: nav
      integer :: msize_

! !REVISION HISTORY:
! 	22May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::msize_'

  msize_=size(nav%displs)
end function msize_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: resize_ - adjust the true size
!
! !DESCRIPTION:
!
!       It is user's responsibility to ensure lsize is no greater than
!   the maximum size of the vector.
!
!       If lsize is not specified, the size of the vector is adjusted
!   to its original size.
!
!
! !INTERFACE:

    subroutine resize_(nav,lsize)
      implicit none
      type(Navigator),intent(inout) :: nav
      integer,optional,intent(in) :: lsize

! !REVISION HISTORY:
! 	22May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::resize_'
  integer :: m

  m=msize_(nav)
  nav%lsize=m
  if(present(lsize)) nav%lsize=lsize
end subroutine resize_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: get_ - get an entry according to the user's preference
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine get_(nav,inav,displ,count,lc,ln,le)
      implicit none
      type(Navigator),intent(in) :: nav
      integer,intent(in) :: inav
      integer,optional,intent(out) :: displ
      integer,optional,intent(out) :: count
      integer,optional,intent(out) :: lc
      integer,optional,intent(out) :: ln
      integer,optional,intent(out) :: le

! !REVISION HISTORY:
! 	22May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::get_'

	! No checking is done here to ensure the entry is a valid.
	! It is the user's responsibility to ensure the index is
	! valid.

  if(present(displ)) displ=nav%displs(inav)
  if(present(count)) count=nav%counts(inav)
  if(present(lc)) lc=nav%displs(inav)+1
  if(present(ln)) ln=nav%counts(inav)
  if(present(le)) le=nav%displs(inav)+nav%counts(inav)

end subroutine get_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_displs_ - referencing to %displs(:)
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_displs_(nav,lbnd,ubnd)
      implicit none
      type(Navigator),intent(in) :: nav
      integer,optional,intent(in) :: lbnd
      integer,optional,intent(in) :: ubnd
      integer,pointer,dimension(:) :: ptr_displs_

! !REVISION HISTORY:
! 	22May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_displs_'
  integer :: lc,le

  if(present(lbnd).or.present(ubnd)) then
    lc=lbound(nav%displs,1)
    if(present(lbnd)) lc=lbnd
    le=ubound(nav%displs,1)
    if(present(ubnd)) le=ubnd
    ptr_displs_ => nav%displs(lc:le)
  else
    le=nav%lsize
    ptr_displs_ => nav%displs(1:le)
  endif

end function ptr_displs_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_counts_ - referencing to %counts(:)
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_counts_(nav,lbnd,ubnd)
      implicit none
      type(Navigator),intent(in) :: nav
      integer,optional,intent(in) :: lbnd
      integer,optional,intent(in) :: ubnd
      integer,pointer,dimension(:) :: ptr_counts_

! !REVISION HISTORY:
! 	22May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_counts_'
  integer :: lc,le

  if(present(lbnd).or.present(ubnd)) then
    lc=lbound(nav%counts,1)
    if(present(lbnd)) lc=lbnd
    le=ubound(nav%counts,1)
    if(present(ubnd)) le=ubnd
    ptr_counts_ => nav%counts(lc:le)
  else
    le=nav%lsize
    ptr_counts_ => nav%counts(1:le)
  endif

end function ptr_counts_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: allgather_ - create a navigator of the all-gathered vector
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine allgather_(localNav,displ,allNav,collNav,comm)
      use m_Collector,only : Collector
      use m_Collector,only : globalSize
      use m_Collector,only :  localSize
      use m_Collector,only : get
      use m_CollectorComm,only : allgatherv
      use m_die,only : die
      implicit none

      type(Navigator),intent(in ) :: localNav	! local Navigator
      integer        ,intent(in ) :: displ	! global displacement
      type(Navigator),intent(out) ::   allNav	! a collected localNav
      type(Collector),intent(in ) :: collNav	! A collector
      integer,intent(in) :: comm		! communicator

! !REVISION HISTORY:
! 	22Aug00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::allgather_'

  integer :: nall,nloc
  integer :: ier

  nloc= localSize(collNav)
	if(nloc/=lsize(localNav)) call die(myname_,	&
		'nloc',nloc,'lsize(localNav)',lsize(localNav))

	! Create a Navigator to _all_ local Navigators

  nall=globalSize(collNav)
  call init_(allNav,nall,stat=ier)
	if(ier/=0) call die(myname_,'init_()',ier)

	! Send local %counts(:) to the Collector

  call allgatherv(localNav%counts      ,allNav%counts,collNav,comm,ier)
	if(ier/=0) call die(myname_,'allgatherv(%counts)',ier)

	! Bump the local %displs(:) by the local vector's global
	! displacement, then send it to the Collector.

  call allgatherv(localNav%displs+displ,allNav%displs,collNav,comm,ier)
	if(ier/=0) call die(myname_,'allgatherv(%displs)',ier)

end subroutine allgather_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: new_ - create an object as a pointer
!
! !DESCRIPTION:
!
! !INTERFACE:

    function new_(mold,stat)
      use m_die, only : die,perr
      use m_mall,only : mall_ison,mall_ci
      implicit none
      type(Navigator),pointer :: mold
      integer,optional,intent(out) :: stat
      type(Navigator),pointer :: new_

! !REVISION HISTORY:
! 	05Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::new_'
  type(Navigator),pointer :: obj
  integer :: ier

  if(present(stat)) stat=0
  allocate(obj,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'allocate()',ier)
	  if(.not.present(stat)) call die(myname_)
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
! !IROUTINE: delete_ - delete an object as a pointer
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine delete_(obj,stat)
      use m_die, only : die,perr
      use m_mall,only : mall_ison,mall_co
      implicit none
      type(Navigator),pointer :: obj
      integer,optional,intent(out) :: stat


! !REVISION HISTORY:
! 	05Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::delete_'
  integer :: ier

  if(present(stat)) stat=0

	if(mall_ison()) call mall_co(1,myname)

  deallocate(obj,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'deallocate()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine delete_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: weighted_define_ - redefine counts by weights
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine weighted_define_(nav,weights,counts)
      use m_die,only : die
      implicit none
      type(Navigator),intent(inout) :: nav
      integer,dimension(:),intent(in) :: weights
      integer,optional,dimension(:),intent(in) :: counts

! !REVISION HISTORY:
! 	12Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::weighted_define_'

  integer :: lc,le,i

	! If redefine counts, do it now.

  if(present(counts)) call define_(nav,counts=counts)

	! redefine counts with weights.  It is an error to call this
	! subroutine without specify counts(:), unless nav has been
	! previously defined.

  do i=1,lsize_(nav)
    call get_(nav,i,lc=lc,le=le)
    nav%counts(i)=sum(weights(lc:le))
  end do

	! redefine displs

  call define_(nav)

end subroutine weighted_define_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: inquire_ - inquire the configuration of an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine inquire_(obj,memsize,ordered,packed)
      use m_die ,only : die
      use m_mall,only : mall_ison,mall_mci,mall_mco
      use m_SortingTools,only : IndexSet
      use m_SortingTools,only : IndexSort
      implicit none
      type(Navigator),intent(in) :: obj
      integer,optional,intent(out) :: memsize	! the upper bound
      logical,optional,intent(out) :: ordered	! in sequence
      logical,optional,intent(out) :: packed	! no "hole"

! !REVISION HISTORY:
! 	13Sep00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::inquire_'
  integer :: memsize_
  logical :: ordered_
  integer :: nnav
  integer :: im,ip,lm,lp
  integer :: i,lc,le,ier
  integer,allocatable,dimension(:) :: indx
  integer,pointer,dimension(:) :: counts
  integer,pointer,dimension(:) :: displs

	! Check the memsize and the order

  memsize_=0
  ordered_=.true.
  nnav=lsize_(obj)

  do i=1,nnav
    call get_(obj,i,lc=lc,le=le)

    if(ordered_   ) ordered_= lc>memsize_
    if(memsize_<le) memsize_= le
  end do

  if(present(ordered)) ordered=ordered_
  if(present(memsize)) memsize=memsize_

	! Check if the navigator is "packed".

  if(present(packed)) then
    packed=.true.
    if(nnav<=0) return

	!________________________________________

	allocate(indx(nnav),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) call mall_mci(indx,myname)
	counts => ptr_counts_(obj)
	displs => ptr_displs_(obj)

    call IndexSet(indx)
    call IndexSort(indx,displs(:),descend=.false.)

	! From 1:memsize, there should be no "hole" if this navigator
	! is "packed".

    lp=indx(1)
    packed= displs(lp)==0

    do im=1,nnav-1	! im and ip are in the sorted orders
      ip=im+1
      lm=indx(im)	! lm and lp are where the two data are.
      lp=indx(ip)

      if(packed) packed= displs(lm)+counts(lm) == displs(lp)
    end do
	!________________________________________

	nullify(displs)
	nullify(counts)
		if(mall_ison()) call mall_mco(indx,myname)
	deallocate(indx,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

  endif

end subroutine inquire_
end module m_Navigator
