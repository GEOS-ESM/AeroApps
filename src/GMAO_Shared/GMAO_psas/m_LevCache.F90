!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_LevCache - Cached gridded data with multiple keys
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_LevCache
      implicit none
      private	! except

      public :: LevCache_init
      public :: LevCache_clean
      public :: LevCache_locate

      interface LevCache_init; module procedure init_; end interface
      interface LevCache_clean; module procedure clean_; end interface
      interface LevCache_locate; module procedure	&
	locate_,	&
	locate2_
      end interface

! !REVISION HISTORY:
! 	05Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_LevCache'

  integer,parameter :: LevCache_MXKEY=2

  integer,save :: LevCache_nkey			! the size of LevCache
  integer,save :: LEvCache_top			! the top of LevCache

  integer,save,allocatable,dimension(:)     :: LevCache_key
  real   ,save,target,allocatable,dimension(:,:,:) :: LevCache_buf

  logical,save :: LevCache_defined=.false.	! the status of LevCache
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize the module
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(nx,ny,nVar)
      use m_mall,only : mall_ison,mall_mci
      use m_die, only : die
      implicit none
      integer,intent(in) :: nx
      integer,intent(in) :: ny
      integer,optional,intent(in) :: nVar

! !REVISION HISTORY:
! 	05Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'

  integer :: ier

  if(LevCache_defined) call die(myname_,'multiple definition')
  LevCache_defined=.true.

  if(ny <=0) call die(myname_,'invalid ny',ny)
  if(nx <=0) call die(myname_,'invalid nx',nx)

  LevCache_nkey=LevCache_MXKEY
  if(present(nVar)) LevCache_nkey=nVar*LevCache_MXKEY

  if(LevCache_nkey <=0) call die(myname_,'invalid nVar',nVar)

  allocate( LevCache_buf(nx,ny,0:LevCache_nkey-1),		&
	    LevCache_key(      0:LevCache_nkey-1), stat=ier	)
	if(ier/=0) call die(myname_,'allocate()',ier)

	if(mall_ison()) then
	  call mall_mci(LevCache_buf,myname)
	  call mall_mci(LevCache_key,myname)
	endif

  LevCache_key(:)=-1
  LevCache_top=0

end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: locate_ - locate a key value in LevCache
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine locate_(key,found,ptr)
      use m_die, only : die
      implicit none
      integer,intent(in) :: key		! the value to be found
      logical,intent(out) :: found	! if the key value is found
      real,pointer,dimension(:,:) :: ptr

! !REVISION HISTORY:
! 	05Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::locate_'
  integer :: itop,ix,i

  if(.not.LevCache_defined) call die(myname_,'undefined')

	! Search the LevCache_key for the given key

  itop=LevCache_top
  ix=-1
  do i=0,LevCache_nkey-1

		! itop isin [0,LevCache_nkey-1]

    if(key==LevCache_key(itop)) then
      ix=itop
      exit
    endif

    itop=mod(itop+1,LevCache_nkey)
  end do

  found=.true.

	! Allocate the next available slot if key is not found.  Note
	! that the cached slots are removed in the order of FIFO, if it
	! is necessary.

  if(ix==-1) then
    found=.false.

    LevCache_top=mod(LevCache_top+1,LevCache_nkey)
    ix=LevCache_top
    LevCache_key(ix)=key
  endif

  ptr=>LevCache_buf(:,:,ix)

end subroutine locate_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: locate2_ - locate two key values in LevCache
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine locate2_(ukey,vkey,found,uptr,vptr)
      use m_die, only : die
      implicit none
      integer,intent(in) :: ukey	! one value to be found
      integer,intent(in) :: vkey	! another value to be found
      logical,intent(out) :: found	! if the key values are found
      real,pointer,dimension(:,:) :: uptr
      real,pointer,dimension(:,:) :: vptr

! !REVISION HISTORY:
! 	05Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::locate2_'
  integer :: itop,iu,iv,i

  if(.not.LevCache_defined) call die(myname_,'undefined')

	! Search the LevCache_key for the given key

  itop=LevCache_top
  iu=-1
  iv=-1
  do i=0,LevCache_nkey-1

		! itop isin [0,LevCache_nkey-1]

    if(ukey==LevCache_key(itop)) iu=itop
    if(vkey==LevCache_key(itop)) iv=itop
    if(iu/=-1 .and. iv/=-1) exit

    itop=mod(itop+1,LevCache_nkey)
  end do

  found=.true.

	! Allocate the next available slot if key is not found.  Note
	! that the cached slots are removed in the order of FIFO, if it
	! is necessary.

  if(iu==-1.or.iv==-1) then
    found=.false.

    if(iu==-1) then
      LevCache_top=mod(LevCache_top+1,LevCache_nkey)
      iu=LevCache_top
      LevCache_key(iu)=ukey
    endif

    if(iv==-1) then
      LevCache_top=mod(LevCache_top+1,LevCache_nkey)
      iv=LevCache_top
      LevCache_key(iv)=vkey
    endif
  endif

  uptr=>LevCache_buf(:,:,iu)
  vptr=>LevCache_buf(:,:,iv)

end subroutine locate2_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean LevCache
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_()
      use m_mall,only : mall_ison,mall_mco
      use m_die, only : die
      implicit none

! !REVISION HISTORY:
! 	05Nov99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

  if(.not.LevCache_defined) call die(myname_,'undefined')
  LevCache_defined=.false.

	if(mall_ison()) then
	  call mall_mco(LevCache_key,myname)
	  call mall_mco(LevCache_buf,myname)
	endif

  deallocate(LevCache_key,LevCache_buf,stat=ier)
	if(ier/=0) call die(myname_,'deallocate()',ier)

  LevCache_nkey=0
  LevCache_top =0

end subroutine clean_
end module m_LevCache
