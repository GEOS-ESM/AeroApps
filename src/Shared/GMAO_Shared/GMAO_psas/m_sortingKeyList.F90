!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_sortingKeyList - a link-list of sortingKeys
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_sortingKeyList
      use m_sortingKeys,only : sortingKeys
      implicit none
      private	! except

      public :: sortingKeyList		! The class data structure
      public :: sortingKeyList_push
      public :: sortingKeyList_clean
      public :: ptr_keySet
      public :: ptr_next

    type sortingKeyList
      type(sortingKeys)   ,pointer,dimension(:) :: keySet
      type(sortingKeyList),pointer :: next
    end type sortingKeyList

    interface sortingKeyList_push ; module procedure	&
	push_  ; end interface
    interface sortingKeyList_clean; module procedure	&
	rclean_; end interface

    interface ptr_keySet; module procedure ptr_keySet_; end interface
    interface ptr_next  ; module procedure ptr_next_  ; end interface

! !REVISION HISTORY:
! 	10Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_sortingKeyList'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: push_ - initialize a sortingKeyList
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine push_(keySet,pList)
      use m_sortingKeys,only : sortingKeys
      use m_mall,only : mall_ison,mall_ci
      use m_die ,only : die
      implicit none
      type(sortingKeys),dimension(:),intent(in) :: keySet
      type(sortingKeyList),pointer :: pList

! !REVISION HISTORY:
! 	10Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::push_'
  type(sortingKeyList),pointer :: next
  integer :: ier
  integer :: lsize

  next => pList
  allocate(pList,stat=ier)
	if(ier/=0) call die(myname_,'allocate(pList)',ier)
	if(mall_ison()) call mall_ci(1,myname)

  pList%next => next

  lsize=size(keySet)
  allocate(pList%keySet(lsize),stat=ier)
	if(ier/=0) call die(myname_,'allocate(%keySet)',ier)
	if(mall_ison()) call mall_ci(lsize,myname)

  pList%keySet=keySet

  nullify(next)
end subroutine push_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rclean_ - recursively clean a sortingKeyList
!
! !DESCRIPTION:
!
! !INTERFACE:

    recursive subroutine rclean_(pList,stat)
      use m_die,only : perr
      use m_mall,only : mall_ison,mall_co
      implicit none
      type(sortingKeyList),pointer :: pList
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	10Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rclean_'
  integer :: ier
  integer :: lsize

  if(present(stat)) stat=0
  if(.not.associated(pList)) return

  call rclean_(pList%next,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'rclean()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

	if(mall_ison()) then
	  lsize=size(pList%keySet)
	  call mall_co(lsize,myname)
	endif
  deallocate(pList%keySet,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'deallocate(%keySet)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

	if(mall_ison()) call mall_co(1,myname)
  deallocate(pList,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'deallocate(pList)',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine rclean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_keySet_ - alias to %keySet
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_keySet_(keyList)
      use m_sortingKeys,only : sortingKeys
      implicit none
      type(sortingKeyList),intent(in) :: keyList
      type(sortingKeys),pointer,dimension(:) :: ptr_keySet_

! !REVISION HISTORY:
! 	10Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_keySet_'

  ptr_keySet_ => keyList%keySet

end function ptr_keySet_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_next_ - alias to pLink%next
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_next_(keyList)
      implicit none
      type(sortingKeyList),intent(in) :: keyList
      type(sortingKeyList),pointer :: ptr_next_

! !REVISION HISTORY:
! 	10Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_next_'

  ptr_next_ => keyList%next

end function ptr_next_

end module m_sortingKeyList
