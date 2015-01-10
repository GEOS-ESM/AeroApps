!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_keyChain - A linked-list of keys
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_keyChain
      implicit none
      private	! except

      public :: keyChain		! The class data structure
      public :: keyChain_init,init	! initialize an object
      public :: clean			! clean an object
      public :: key			! the current key value
      public :: ptr_next		! the next keyChain level

    type keyChain
      private
      integer :: key
      type(keyChain),pointer :: next
    end type keyChain

    interface keyChain_init; module procedure rinit_; end interface
    interface init ; module procedure rinit_ ; end interface
    interface clean; module procedure rclean_; end interface

    interface key  ; module procedure key_  ; end interface
    interface ptr_next; module procedure ptr_next_; end interface

! !REVISION HISTORY:
! 	24Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_keyChain'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rinit_ - initialize an object, recursively
!
! !DESCRIPTION:
!
! !INTERFACE:

    recursive subroutine rinit_(keyCh,keys)
      use m_die ,only : die
      use m_mall,only : mall_ison,mall_ci
      implicit none
      type(keyChain)      ,intent(out) :: keyCh	! keys in a linked-list
      integer,dimension(:),intent(in)  :: keys	! one or more keys

! !REVISION HISTORY:
! 	24Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rinit_'
  integer :: ier

	! At least, there should be one key, since keyCh is not a
	! pointer, but an actual variable with keyCh%key accessible.

  if(size(keys)<=0) call die(myname_,'invalid size(keys)',size(keys))

  keyCh%key=keys(1)	! set this node
  nullify(keyCh%next)

  if(size(keys)>1) then

    allocate(keyCh%next,stat=ier)
	if(ier/=0) call die(myname_,'allocate()',ier)
	if(mall_ison()) call mall_ci(1,myname)

    call rinit_(keyCh%next,keys(2:))
  endif

end subroutine rinit_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rclean_ - clean an object, recursively
!
! !DESCRIPTION:
!
! !INTERFACE:

    recursive subroutine rclean_(keyCh)
      use m_die ,only : die
      use m_mall,only : mall_ison,mall_co
      implicit none
      type(keyChain),intent(inout) :: keyCh

! !REVISION HISTORY:
! 	24Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rclean_'
  integer :: ier

  if(associated(keyCh%next)) then
    call rclean_(keyCh%next)
	if(mall_ison()) call mall_co(1,myname)
    deallocate(keyCh%Next,stat=ier)
	if(ier/=0) call die(myname_,'deallocate()',ier)
  endif

end subroutine rclean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: key_ - get the key value
!
! !DESCRIPTION:
!
! !INTERFACE:

    function key_(keyCh)
      implicit none
      type(keyChain),intent(in) :: keyCh
      integer :: key_

! !REVISION HISTORY:
! 	24Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::key_'

  key_=keyCh%key

end function key_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_next_ - alias to keyCh%next
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_next_(keyCh)
      implicit none
      type(keyChain),intent(in) :: keyCh
      type(keyChain),pointer :: ptr_next_

! !REVISION HISTORY:
! 	24Oct00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_next_'

  ptr_next_ => keyCh%next

end function ptr_next_
end module m_keyChain
